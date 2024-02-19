.. _API_2_sim:

***************
steps.API_2.sim
***************

The ``sim`` module contains classes related to both serial and parallel simulations.

================
Simulation paths
================

Values can be set and retrieved from the simulation by using `SimPaths <API_sim.rst#steps.API_2.sim.SimPath>`_ from the :py:class:`Simulation` object to the desired value.
The following tool allows users to know which paths can be set or accessed:

.. raw:: html

    <div>
        <div style="display:inline-block">Get / Set</br><select style="width:100px" class="SimPathSelect" id="GetSet"></select></div>
        <div style="display:inline-block">Solver</br><select style="width:200px" class="SimPathSelect" id="Solver"></select></div>
    </div>
    <div>
        <div style="display:inline-block">Location</br><select style="width:200px;height:200px;" class="SimPathSelect" id="Location" multiple></select></div>
        <div style="display:inline-block">Object / Property</br><select style="width:200px;height:200px;" class="SimPathSelect" id="Item" multiple></select></div>
        <div style="display:inline-block">Property</br><select style="width:200px;height:200px;" class="SimPathSelect" id="Property" multiple></select></div>
        <h4>Examples:</h4>
        <div class="DisplayBox line-block" id="SimPath" style="font-size:25px;line-height:30px;width:610px;min-height:100px;"></div>
        <h4>Description:</h4>
        <div class="DocDisplay line-block" id="SimPathDoc" style="width:610px;min-height:200px;color:#666666;"></div>
    </div>

======================
Detailed documentation
======================

* Simulation

 * :py:class:`steps.API_2.sim.Simulation`
 * :py:class:`steps.API_2.sim.SBMLSimulation`

* Simulation utilities

 * :py:class:`steps.API_2.sim.SimPath`
 * :py:class:`steps.API_2.sim.VesiclePathReference`
 * :py:class:`steps.API_2.sim.VesicleReference`
 * :py:class:`steps.API_2.sim.RaftReference`
 * :py:class:`steps.API_2.sim.VesicleList`
 * :py:class:`steps.API_2.sim.RaftList`
 * :py:class:`steps.API_2.sim.LinksSpecList`

* Parallel simulations

 * :py:class:`steps.API_2.sim.MPI`

----------

.. automodule:: steps.API_2.sim
    :members:
    :undoc-members:
    :special-members:
    :show-inheritance:



