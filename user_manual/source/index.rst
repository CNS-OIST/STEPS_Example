
.. raw:: html

    <img src="_static/logo.svg" style="background:#2980b9;padding:10px;border-radius:5px"/>

----------

STEPS is a package for exact stochastic simulation of reaction-diffusion systems in arbitrarily complex 3D geometries. Our core simulation algorithm is an implementation of Gillespie's SSA, extended to deal with diffusion of molecules over the elements of a 3D tetrahedral mesh.


While it was mainly developed for simulating detailed models of neuronal signaling pathways in dendrites and around synapses, it is a general tool and can be used for studying any biochemical pathway in which spatial gradients and morphology are thought to play a role.


STEPS also supports accurate and efficient computation of local membrane potentials on tetrahedral meshes, with the addition of voltage-gated channels and currents. Tight integration between the reaction-diffusion calculations and the tetrahedral mesh potentials allows detailed coupling between molecular activity and local electrical excitability.

We have implemented STEPS as a set of Python modules, which means STEPS users can use Python scripts to control all aspects of setting up the model, generating a mesh, controlling the simulation and generating and analyzing output. The core computational routines are still implemented as C/C++ extension modules for maximal speed of execution.

STEPS 3.0.0 and above provide early parallel solution for stochastic spatial reaction-diffusion and electric field simulation.

STEPS 3.6.0 and above provide a new set of APIs (API2) to speedup STEPS model development. Models developed with the old API (API1) are still supported.

STEPS 4.0.0 and above provide a distributed solution for stochastic spatial reaction-diffusion and electric field simulation.

----------

STEPS User Manual and API References
====================================

Contents:

.. toctree::
   :maxdepth: 2
   :numbered:
  
   getting_started.ipynb
   Python_interfaces.ipynb
   API_2/Interface_Tutorial_1_wm.ipynb
   API_2/Interface_Tutorial_2_IP3.ipynb
   API_2/Interface_Tutorial_3_Diffusion.ipynb
   API_2/Interface_Tutorial_3.1_Paraview.ipynb
   API_2/Interface_Tutorial_4_Complexes.ipynb
   API_2/Interface_Tutorial_5_Efield.ipynb
   API_2/Interface_Tutorial_6_MPI.ipynb
   API_2/Interface_Tutorial_7_visual.ipynb
   API_2/Interface_Tutorial_8_Caburst.ipynb
   API_1_guide.rst
   API_2/API_ref.rst
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

