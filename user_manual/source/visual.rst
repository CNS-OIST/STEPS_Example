.. _visual:

****************************************
Visualization Toolkit
****************************************

The simulation scripts described in this chapter are available at `STEPS_Example repository <https://github.com/CNS-OIST/STEPS_Example/tree/master/publication_models/Chen_FNeuroinf_2014>`_.

In this chapter, we'll use simple models as examples to introduce the use of visualization toolkit described in `Python-based geometry preparation and simulation visualization toolkits for STEPS <http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00037/abstract>`_.

API reference of this module can be accessed via :ref:`API_visual`

Prerequisites
===================

The following third-party python packages are required to use the visualization toolkit.

1. `NumPy <http://www.numpy.org/>`_
2. `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_
3. `PyQtGraph <http://www.pyqtgraph.org/>`_

Installation
===================

If the prerequsite packages are installed, the visualization toolkit will be installed along with
STEPS following the standard installation.

The toolkit module can be import in Python by ::

    import steps.visual

Examples
=========

Visualization of Spatial IP3 Receptor Model
-------------------------------------------

Here we use the IP3 receptor model described in :ref:`ip3` and its spatial extension in :ref:`spatial_ip3`
to demostrate how to use the visualization toolkit in practice.

Biochemical Model
^^^^^^^^^^^^^^^^^

We stored the biochemical model as a individual python script, **ip3r_model.py**, ::

    # IP3 receptor model

    import steps.model as smodel
    import steps.geom as sgeom

    ###############################################################################

    def getModel():
        mdl = smodel.Model()
        
        # chemical species objects
        Ca = smodel.Spec('Ca', mdl)				# Calcium
        IP3 = smodel.Spec('IP3', mdl)			# IP3
        
        # receptor state objects
        R = smodel.Spec('R', mdl)				# IP3 receptor in 'naive' state
        RIP3 = smodel.Spec('RIP3', mdl)			# bound IP3 
        Ropen = smodel.Spec('Ropen', mdl)		# bound IP3 and Ca (open)
        RCa = smodel.Spec('RCa', mdl)			# 1 bound Ca to inactivation site
        R2Ca = smodel.Spec('R2Ca', mdl)			# 2 bound Ca to inactivation sites
        R3Ca = smodel.Spec('R3Ca', mdl)			# 3 bound Ca to inactivation sites
        R4Ca = smodel.Spec('R4Ca', mdl)			# 4 bound Ca to inactivation sites
        
        surfsys = smodel.Surfsys('ssys', mdl)
        
        # The 'forward' binding reactions: 
        R_bind_IP3_f = smodel.SReac('R_bind_IP3_f', surfsys, \
            olhs=[IP3], slhs=[R], srhs=[RIP3])
        RIP3_bind_Ca_f = smodel.SReac('RIP3_bind_Ca_f', surfsys, \
            olhs=[Ca], slhs=[RIP3], srhs = [Ropen])
        R_bind_Ca_f = smodel.SReac('R_bind_Ca_f', surfsys, \
            olhs=[Ca], slhs=[R], srhs=[RCa])
        RCa_bind_Ca_f = smodel.SReac('RCa_bind_Ca_f', surfsys, \
            olhs=[Ca], slhs=[RCa],srhs = [R2Ca])
        R2Ca_bind_Ca_f = smodel.SReac('R2Ca_bind_Ca_f', surfsys, \
            olhs=[Ca], slhs= [R2Ca], srhs = [R3Ca])
        R3Ca_bind_Ca_f = smodel.SReac('R3Ca_bind_ca_f', surfsys, \
            olhs=[Ca], slhs=[R3Ca], srhs=[R4Ca])
            
        # The 'backward' binding reactions:
        R_bind_IP3_b = smodel.SReac('R_bind_IP3_b', surfsys, \
            slhs=[RIP3], orhs=[IP3], srhs=[R])
        RIP3_bind_Ca_b = smodel.SReac('RIP3_bind_Ca_b', surfsys, \
            slhs=[Ropen], orhs=[Ca], srhs=[RIP3])
        R_bind_Ca_b = smodel.SReac('R_bind_Ca_b', surfsys, \
            slhs=[RCa], orhs=[Ca], srhs=[R])
        RCa_bind_Ca_b = smodel.SReac('RCa_bind_Ca_b', surfsys, \
            slhs=[R2Ca], orhs=[Ca], srhs=[RCa])
        R2Ca_bind_Ca_b = smodel.SReac('R2Ca_bind_Ca_b', surfsys, \
            slhs=[R3Ca], orhs=[Ca], srhs= [R2Ca])
        R3Ca_bind_Ca_b = smodel.SReac('R3Ca_bind_ca_b', surfsys, \
            slhs=[R4Ca], orhs=[Ca], srhs=[R3Ca])
        
        # Ca ions passing through open IP3R channel
        R_Ca_channel_f = smodel.SReac('R_Ca_channel_f', surfsys, \
            ilhs=[Ca], slhs=[Ropen], orhs=[Ca], srhs=[Ropen])
        R_Ca_channel_b = smodel.SReac('R_Ca_channel_b', surfsys, \
            olhs=[Ca], slhs=[Ropen], irhs=[Ca], srhs=[Ropen])
        
        # The reaction constants
        R_bind_IP3_f.setKcst(1000e6)
        R_bind_IP3_b.setKcst(25800)
        RIP3_bind_Ca_f.setKcst(8000e6)
        RIP3_bind_Ca_b.setKcst(2000)
        R_bind_Ca_f.setKcst(8.889e6)
        R_bind_Ca_b.setKcst(5)
        RCa_bind_Ca_f.setKcst(20e6)
        RCa_bind_Ca_b.setKcst(10)
        R2Ca_bind_Ca_f.setKcst(40e6)
        R2Ca_bind_Ca_b.setKcst(15)
        R3Ca_bind_Ca_f.setKcst(60e6)
        R3Ca_bind_Ca_b.setKcst(20)
        
        # Corresponds to Ca input ~ 20000/ms for open receptor
        R_Ca_channel_f.setKcst(8e6)          
        R_Ca_channel_b.setKcst(8e6)           
        
        return mdl

Typical STEPS Simulation Routine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the simulation script, we first import the biocemical model and some standard STEPS modules.
We also define the diffuson constants according to publication. ::

    # IP3 receptor mesh simulation

    import steps.model as smodel
    import steps.geom as swm
    import steps.rng as srng
    import steps.solver as ssolver

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    # DIFFUSION

    # Source:
    #   Allbritton, N.L., Meyer, T., and Stryer, L. (1992). 
    #   Range of messenger action of calcium ion and inositol 
    #   1,4,5-triphosphate. Science 258, 1812-1815.
    DCST_Ca = 0.065e-9
    DCST_IP3 = 0.283e-9

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    import ip3r_model 

    # Import model
    mdl = ip3r_model.getModel()

We then add a :mod:`steps.model.Volsys` volume system to the model to host the :mod:`steps.model.Diff` 
diffution rules for Calcium and IP3. ::

    volsys = smodel.Volsys('vsys', mdl)

    # Fetch reference to Calcium and IP3 Spec objects
    Ca = mdl.getSpec('Ca')
    IP3 = mdl.getSpec('IP3')

    # Create diffusion rules
    Ca_diff = smodel.Diff('Ca_diff', volsys, Ca, DCST_Ca)
    IP3_diff = smodel.Diff('IP3_diff', volsys, IP3, DCST_IP3)

Now we load the :mod:`steps.geom.Tetmesh` from file prepared in :ref:`spatial_ip3`, using :func:`steps.utilities.meshio.loadMesh`. ::

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Import mesh
    import steps.utilities.meshio as meshio

    mesh = meshio.loadMesh("ip3r_mesh")[0]

We then create the random number generator and the :mod:`steps.solver.Tetexact` solver,
and initialize the simulation by adding molecules into compartments and patch. ::

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    # Create random number generator
    r = srng.create('mt19937', 512)
    r.initialize(456)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    # Create reaction-diffusion solver object
    sim = ssolver.Tetexact(mdl, mesh, r)

    # Setup initial condition
    sim.setCompConc('cyt', 'Ca', 3.30657e-8)
    sim.setCompConc('cyt', 'IP3', 2.5e-6)
    sim.setCompConc('ER', 'Ca', 150e-6)
    sim.setPatchCount('memb', 'R', 16)

The above scripts are typical STEPS simulation routines. With the model, geometry and simulation
solver ready, we can now work on constructing the visualization system.

Visualization Routine
^^^^^^^^^^^^^^^^^^^^^

First, we import the visualization module and pyqtgraph module, also create a standard QApplication instance. ::

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Visualization
    import pyqtgraph as pg
    import steps.visual as visual

    # Visualization initialization
    app = pg.mkQApp()

Plot Display
''''''''''''

Now let's create a plot display so that we can have a quantitative view of the simulation.
We first create a :mod:`steps.visual.PlotDisplay` instance, which will be the host of all our plots.
::

    # Create plot display
    plots = visual.PlotDisplay("IP3 Receptor Model", size = (600, 400))

The :mod:`steps.visual.PlotDisplay` class provides varies functions for different plotting requirements.
For example, :func:`steps.visual.PlotDisplay.addCompSpecPlot` adds a plot to the display, which displays 
molecule count/Concentration changes during simulation. We can also setup different display features 
for the plot, such as axis labels, data style, etc. The following code creates a plot showing 
calcium Concentration changes in cytosol. ::

    # Create Plots
    pen = pg.mkPen(color=(255,255,255), width=2)
    p = plots.addCompSpecPlot("<span style='font-size: 16pt'>Ca_cyt", sim, "cyt", "Ca", data_size = 1000,y_range= [0, 1e-5], measure = "conc", pen=(255, 0.647 * 255, 0))
    p.getAxis('left').setPen(pen)
    p.getAxis('bottom').setPen(pen)
    p.showGrid(x=True, y=True)
    labelStyle = {'color': '#ffffff', 'font-size': '16px'}
    p.setLabel('bottom', 'Time', 's', **labelStyle)

The plot display arrange plot items in a grid system. As we've completed the above plot, we switch to
next row and start creating a new plot which displays molecule changes of "Ropen" species (that is the IP3 
receptor at open state) on the membrane patch. ::

    plots.nextRow()

    p = plots.addPatchSpecPlot("<span style='font-size: 16pt'>Ropen_memb", sim, "memb", "Ropen", data_size = 1000,y_range= [0, 10], pen=(255, 0, 255))
    p.getAxis('left').setPen(pen)
    p.getAxis('bottom').setPen(pen)
    p.showGrid(x=True, y=True)
    p.setLabel('bottom', 'Time', 's', **labelStyle)

Simulation Display
''''''''''''''''''

We start working on the actual simulation displays. In this example, we would like to create multiple
display windows, one for overview of the complete system, and several others for varies components.

A simulation display can be constructed by creating a :mod:`steps.visual.SimDisplay` object. ::

    # Create simulation displays
    full_display = visual.SimDisplay("Full View", w = 600, h = 400)
    ER_display = visual.SimDisplay("ER", w = 600, h = 400)
    cytIP3_display = visual.SimDisplay("Cyt IP3", w = 600, h = 400)
    cytCa_display = visual.SimDisplay("Cyt Calcium", w = 600, h = 400)
    memb_display = visual.SimDisplay("memb", w = 600, h = 400)

Now it is time to add different visual components to the displays.
The Visualization toolkit provides two major types of visual components, static and dynamic.

Static components include :mod:`steps.visual.VisualCompMesh` for Visualizing compartment mesh, and
:mod:`steps.visual.VisualPatchMesh` for Visualizing patch mesh. We create two :mod:`steps.visual.VisualCompMesh` instances for both cytosol and ER, and a :mod:`steps.visual.VisualPatchMesh` instance for ER membrane. ::

    # Create static mesh components
    ER_view = visual.VisualCompMesh("ER", full_display, mesh, "ER", color = [0.678, 1.000, 0.184, 0.05])
    cyt_view = visual.VisualCompMesh("cyt", full_display, mesh, "cyt", color = [0.941, 1.000, 0.941, 0.05])
    memb_view = visual.VisualPatchMesh("memb", full_display, mesh, "memb", color = [1.000, 0.973, 0.863, 0.05])

Dynamic components include several variations of molecule visualization components for species in compartments
or on patches. Here we use two different components, :mod:`steps.visual.VisualCompSpec` for visualizing 
compartmental species such as calcium in ER and cytosol, as well as IP3 in cytosol. And :mod:`steps.visual.VisualPatchChannel` for visualizing different states of IP3 receptors on ER membrane. ::

    # Create dynamic species components
    Ca_ER = visual.VisualCompSpec("Ca_ER", full_display, mesh, sim, "ER", "Ca", [1.000, 0.647, 0.000, 1.0], spec_size = 0.005)
    IP3_cyt = visual.VisualCompSpec("IP3_cyt", full_display, mesh, sim, "cyt", "IP3", [1.0, 0.0, 0.0, 1.0], spec_size = 0.005)
    Ca_cyt = visual.VisualCompSpec("Ca_cyt", full_display, mesh, sim, "cyt", "Ca", [1.000, 0.647, 0.000, 1.0], spec_size = 0.005)
    IP3R_MEMB = visual.VisualPatchChannel("IP3R_memb", full_display, mesh, sim, "memb", {"R" : [0.0, 0.0, 1.0, 1.0], "RIP3" : [1.0, 0.0, 1.0, 0.2], "Ropen" : [1.0, 0.0, 1.0, 1.0], "RCa" : [0.0, 0.0, 1.0, 0.8], "R2Ca" : [0.0, 0.0, 1.0, 0.6], "R3Ca" : [0.0, 0.0, 1.0, 0.4], "R4Ca" : [0.0, 0.0, 1.0, 0.2]}, spec_size = 0.01)

We then add these visual components to associated simulation displays ::

# Add associated components to individual displays

    ER_display.addItem(ER_view)
    ER_display.addItem(Ca_ER)

    cytCa_display.addItem(cyt_view)
    cytCa_display.addItem(Ca_cyt)

    cytIP3_display.addItem(cyt_view)
    cytIP3_display.addItem(IP3_cyt)

    memb_display.addItem(memb_view)
    memb_display.addItem(IP3R_MEMB)
    

Simulation Control
''''''''''''''''''
The final task is to create a :mod:`steps.visual.SimControl` controller and assigned simulation and displays
to it. ::

    # Add simulation and displays to control
    x = visual.SimControl([sim], [ER_display, cytIP3_display, cytCa_display, memb_display, full_display],[plots], end_time= 1.0, upd_interval = 0.0001)

    # Enter visualization loop
    app.exec_()

Showcase: Plots and Visual Components
===========================================

One essential step of STEPS simulation visualization is to choose suitable plotting function and visual component
based on project requirement. To provide an intuitive concept of each component, here we showcase some examples 
of them in practice.

Plots
-----

Count/Concentration Plot
^^^^^^^^^^^^^^^^^^^^^^^^^
    * :mod:`steps.visual.PlotDisplay.addCompSpecPlot`
    * :mod:`steps.visual.PlotDisplay.addTetsSpecPlot`
    * :mod:`steps.visual.PlotDisplay.addPatchSpecPlot`
    * :mod:`steps.visual.PlotDisplay.addTrisSpecPlot`
    * :mod:`steps.visual.PlotDisplay.addCompSumSpecsPlot`
    * :mod:`steps.visual.PlotDisplay.addPatchSumSpecsPlot`
    
.. raw:: html

        <object width="640" height="480"><param name="movie"
        value="http://www.youtube.com/v/5gv7wRIGRhM"></param><param
        name="allowFullScreen" value="true"></param><param
        name="allowscriptaccess" value="always"></param><embed
        src="http://www.youtube.com/v/5gv7wRIGRhM"
        type="application/x-shockwave-flash" allowscriptaccess="always"
        allowfullscreen="true" width="640"
        height="480"></embed></object>
        
Distribution Plot
^^^^^^^^^^^^^^^^^
    * :mod:`steps.visual.PlotDisplay.addCompSpecDist`
    * :mod:`steps.visual.PlotDisplay.addPatchSpecDist`
    * :mod:`steps.visual.PlotDisplay.addTetsSpecDist`
    * :mod:`steps.visual.PlotDisplay.addTrisSpecDist`
    * :mod:`steps.visual.PlotDisplay.addROISpecDist`

.. raw:: html

        <object width="640" height="480"><param name="movie"
        value="http://www.youtube.com/v/sb67Cj8u3A0"></param><param
        name="allowFullScreen" value="true"></param><param
        name="allowscriptaccess" value="always"></param><embed
        src="http://www.youtube.com/v/sb67Cj8u3A0"
        type="application/x-shockwave-flash" allowscriptaccess="always"
        allowfullscreen="true" width="640"
        height="480"></embed></object>

Visual Component
----------------

Static Component
^^^^^^^^^^^^^^^^
    * :mod:`steps.visual.VisualCompMesh`
    * :mod:`steps.visual.VisualPatchMesh`

.. raw:: html

        <object width="640" height="480"><param name="movie"
        value="http://www.youtube.com/v/w457Cv-vJdI"></param><param
        name="allowFullScreen" value="true"></param><param
        name="allowscriptaccess" value="always"></param><embed
        src="http://www.youtube.com/v/w457Cv-vJdI"
        type="application/x-shockwave-flash" allowscriptaccess="always"
        allowfullscreen="true" width="640"
        height="480"></embed></object>

Dynamic Component
^^^^^^^^^^^^^^^^^
Diffusive Species
''''''''''''''''''
    * :mod:`steps.visual.VisualTetsSpec`
    * :mod:`steps.visual.VisualCompSpec`
    * :mod:`steps.visual.VisualROITetsSpec`
    * :mod:`steps.visual.VisualTrisSpec`
    * :mod:`steps.visual.VisualPatchSpec`
    * :mod:`steps.visual.VisualROITrisSpec`

.. raw:: html

        <object width="640" height="480"><param name="movie"
        value="http://www.youtube.com/v/xvnBUJxoU7Y"></param><param
        name="allowFullScreen" value="true"></param><param
        name="allowscriptaccess" value="always"></param><embed
        src="http://www.youtube.com/v/xvnBUJxoU7Y"
        type="application/x-shockwave-flash" allowscriptaccess="always"
        allowfullscreen="true" width="640"
        height="480"></embed></object>

Non-diffusive Channel Species
'''''''''''''''''''''''''''''
    * :mod:`steps.visual.VisualTrisChannel`
    * :mod:`steps.visual.VisualPatchChannel`
    * :mod:`steps.visual.VisualROITrisChannel`
    
.. raw:: html

        <object width="640" height="480"><param name="movie"
        value="http://www.youtube.com/v/Zv_oFUfVJk0"></param><param
        name="allowFullScreen" value="true"></param><param
        name="allowscriptaccess" value="always"></param><embed
        src="http://www.youtube.com/v/Zv_oFUfVJk0"
        type="application/x-shockwave-flash" allowscriptaccess="always"
        allowfullscreen="true" width="640"
        height="480"></embed></object>


