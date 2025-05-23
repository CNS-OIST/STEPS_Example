.. _API_1_ref:

****************
API_1 References
****************

Namespaces:

.. toctree::
   :maxdepth: 3
   
   API_model.rst
   API_geom.rst
   API_rng.rst
   API_solver.rst
   API_mpisolver.rst
   API_utilities_meshio.rst
   API_utilities_collections.rst
   API_utilities_sbml.rst
   API_visual.rst
   
Currently STEPS has 6 major namespaces, divided by functionality. 
The :mod:`steps.API_1.model` namespace defines classes of chemical species in the model, 
as well as the reaction and diffusion objects. The :mod:`steps.API_1.geom` namespace 
contains classes that define the geometry of the model (i.e. tetrahedral meshes, 
compartments and patches). The :mod:`steps.API_1.rng` namespace contains the random number 
generators (currently only one generator is available). The :mod:`steps.API_1.solver` and :mod:`steps.API_1.mpi.solver` namespaces 
implement the simulation solvers. The :mod:`steps.API_1.utilities` namespace provides useful 
tools related to the model and geometry description, currently providing 3D 
tetrahedral mesh loading and saving support. 

.. figure:: ../images/steps.jpg
   :align:  center
   
   STEPS main class diagram.

**Methods and attributes**

For most of the set and get methods, STEPS provides corresponding attributes 
which can be accessed and modified directly. For example (spec is a reference 
to an object of type :class:`steps.API_1.model.Spec`)::

    >>> spec.setID('Ca')
    
is equivalent to::

    >>> spec.id = 'Ca'
    
and,

::

    >>> spec_id = spec.getID()
    
is equivalent to::

    >>> spec_id = spec.id
    
Available methods and attributes of an object can be found by the built-in 
dir command in Python. A convention maintained in STEPS is that all attributes 
are completely lower-case::

    >>> dir(spec)
    ['__class__', '__del__', '__delattr__', '__dict__', '__doc__', '__getattr__', 
    '__getattribute__', '__hash__', '__init__', '__module__', '__new__',
    '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__str__', 
    '__swig_destroy__', '__swig_getmethods__', '__swig_setmethods__', 
    '__weakref__', 'getID', 'getModel', 'id', 'model', 'setID', 'this']
    
All attributes beginning with double underscore are built-in attributes related 
to the implementation of the type and can be ignored. So, for our 
:func:`steps.API_1.model.Spec` object (one of the most basic objects in STEPS) we have 
functions ``getID``, ``setID`` and ``getModel`` available with attributes ``id`` 
and ``model`` related to those functions.


