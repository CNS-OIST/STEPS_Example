.. _API_2_sim:

***************
steps.API_2.sim
***************

* Simulation

 * :py:class:`steps.API_2.sim.Simulation`
 * :py:class:`steps.API_2.sim.SBMLSimulation`

* Simulation utilities

 * :py:class:`steps.API_2.sim.SimPath`

* Parallel simulations

 * :py:class:`steps.API_2.sim.MPI`


================
Simulation paths
================

Values can be set and retrieved from the simulation by using `SimPaths <API_sim.rst#steps.API_2.sim.SimPath>`_ from the simulation object to the desired value.
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
        <code class="DisplayBox" id="SimPath" style="font-size:25px;width:610px;"></code>
        <h4>Description:</h4>
        <div class="DocDisplay" id="SimPathDoc" style="width:610px;height:200px;color:#666666;"></div>
    </div>

----------

.. automodule:: steps.API_2.sim
    :members:
    :undoc-members:
    :inherited-members: ndarray
    :special-members:



