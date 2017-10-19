.. _preliminary_func:

***************************
Preliminary Functionalities
***************************

STEPS has implemented several functionalities which are still preliminary.
Most of them are written in Python scripts so that user can extend them easily
to meet their individual requirements.

.. Note:: 

    We welcome any modification and/or extension of these functionarities. 
    If you would like to share your extensions with other users, 
    please contact us steps.dev@gmail.com

meshio: Mesh Input/Output Module
================================

Tetrahedral mesh based simulation is one of the important features in STEPS.
As there are many tetrahedral mesh generators available in public, STEPS provides
a python based meshio module for importing meshes from 3rd party mesh files.

Currently meshio module also provides one-stop importing functions for Abaqus and Tetgen
formats.

To import mesh from Abaqus(\*.inp) files::
    
    mesh = steps.utilities.meshio.importAbaqus(pathroot, scale)
    
To import mesh from TetGen files::

    mesh = steps.utilities.meshio.importTetGen(pathroot)

STEPS also provides methods to save any Tetrahedral mesh described in a :class:`steps.geom.Tetmesh`
object in a STEPS format, which contains a lot of connectivity information that is absent from other formats
yet is vital for a STEPS simulation. Once a mesh has been imported once from an outside package the mesh may 
be saved in this format and, because all connectivity information is saved, mesh loading time is reduced 
considerably. Therefore it is recommended to save all large meshes in STEPS format after importing from
the mesh-generator the first time. 

To save a mesh into STEPS format::

   steps.utilities.meshio.saveMesh(pathroot, tetmesh)

To load a mesh from STEPS format::

    mesh, comps, patches = steps.utilities.meshio.loadMesh(pathroot)
    
To extend this module, experienced user could modify ``steps/utilities/meshio.py`` 
in the installed package.

More details of the meshio module can be found in :doc:`API_utilities_meshio`.

Checkpointing
=============

Currently STEPS provides a basic checkpointing functionality
for :class:`steps.solver.Wmdirect` and :class:`steps.solver.Tetexact` via 
Python's CPickle module. 

.. note::
    
    Currently STEPS only checkpoints the species distributions in the simulation,
    this means the user will need to reset any clamped or deactived reaction/diffusion
    when restoring from checkpointing files.

To run a simulation with automatic checkpointing, simply add the cp_interval as
the second parameter in ``run`` or ``advance`` function::

    sim.run(endtime, cp_interval)
    sim.advance(adv, cp_interval)

Experienced user can extend this functionality by modifying ``steps/solver.py``
in the installed package.


