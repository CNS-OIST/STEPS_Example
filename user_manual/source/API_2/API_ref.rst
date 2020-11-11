.. _API_2_ref:

****************
API_2 References
****************

Using API_2
===========

In order to use API_2, python scripts must start with the ``import steps.interface`` directive. 
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
