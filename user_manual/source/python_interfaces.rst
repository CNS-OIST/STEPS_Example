.. _interfaces:

**************************
Python interfaces to STEPS
**************************

Since version 3.6, STEPS is shipped with a new, more pythonic interface. In this documentation, we will refer to the classical interface as API_1 and to the new one as API_2.

Compatibility between APIs
==========================

Users have to decide whether they want to use API_1 or API_2, it is not advised to try to mix both in a single python script.
Models that were written using API_1 are still usable without any changes.
Users that want to use API_2 have to explicitely specify it at the beginning of their python scripts (see :doc:`API_2/API_ref`).

Guides
======

Chapters 3 to 12 cover the use of the classical interface and introduce STEPS concepts in details along the way:

    - :doc:`well_mixed`
    - :doc:`ip3`
    - :doc:`diffusion`
    - :doc:`surface_diffusion`
    - :doc:`diffusion_boundary`
    - :doc:`surface_diffusion_boundary`
    - :doc:`memb_pot`
    - :doc:`stoch_spikes`
    - :doc:`parallel`
    - :doc:`visual`

Users that are already familiar with STEPS concepts and the classical interface can learn how to use the new interface in `chapter 13 <API_2/guide.rst>`_:

    - :doc:`API_2/Interface_Tutorial_1_wm`
    - :doc:`API_2/Interface_Tutorial_2_IP3`
    - :doc:`API_2/Interface_Tutorial_3_Diffusion`
    - :doc:`API_2/Interface_Tutorial_4_Complexes`
    - :doc:`API_2/Interface_Tutorial_5_Efield`
    - :doc:`API_2/Interface_Tutorial_6_MPI`
    - :doc:`API_2/Interface_Tutorial_7_visual`

API references
==============

Detailed documentations are available for both APIs:
    - :doc:`API_1/API_ref`
    - :doc:`API_2/API_ref`
