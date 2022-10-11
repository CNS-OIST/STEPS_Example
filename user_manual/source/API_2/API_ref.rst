.. _API_2_ref:

**************
API References
**************

In order to use the most recent STEPS API, python scripts must start with the ``import steps.interface`` directive. 
STEPS modules should preferentially be imported using the ``from steps.X import *`` syntax but no 
other modules should be imported in this way in order to avoid name collisions.

::

    import steps.interface

    from steps.model import *
    from steps.geom import *
    from steps.rng import *
    from steps.sim import *
    from steps.saving import *
    from steps.utils import *
    from steps.visual import *

Note that some functionalites of STEPS API, like `auto-naming <API_utils.rst#steps.API_2.utils.NamedObject.Create>`_, will not work in an interactive python shell but will work as expected in a python script or a jupyter notebook.

The currently used API can be get / set with functions from the ``steps`` module:

.. automodule:: steps
    :members:

----------

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   API_model
   API_geom
   API_rng
   API_saving
   API_sim
   API_utils
   API_visual
   ../API_1/API_ref
