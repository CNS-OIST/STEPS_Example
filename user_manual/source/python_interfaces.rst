.. _interfaces:

**************************
Python interfaces to STEPS
**************************

Since version 3.6, STEPS is shipped with a new, more pythonic interface. In this documentation, we will refer to the original interface as API_1 and to the new one as API_2. Unless mentioned explicitely, we will use API_2 by default.

Compatibility between APIs
==========================

Models that were written using API_1 are still usable without any changes. Users are however encouraged to start using API_2 for new models. It is not advised to try to mix both APIs in a single python script.
Until API_1 becomes deprecated, users that want to use API_2 have to explicitely specify it at the beginning of their python scripts (see :doc:`API_2/API_ref`).

Guides
======

Chapters 3 to 9 introduce STEPS concepts in details through examples:

    - :doc:`API_2/Interface_Tutorial_1_wm`
    - :doc:`API_2/Interface_Tutorial_2_IP3`
    - :doc:`API_2/Interface_Tutorial_3_Diffusion`
    - :doc:`API_2/Interface_Tutorial_4_Complexes`
    - :doc:`API_2/Interface_Tutorial_5_Efield`
    - :doc:`API_2/Interface_Tutorial_6_MPI`
    - :doc:`API_2/Interface_Tutorial_7_visual`
    - :doc:`API_2/Interface_Tutorial_8_CaBurst`

Equivalent chapters using API_1 are still available for reference (see :doc:`API_1_guide`).

Finally, :doc:`API_2/API_ref` provides detailed documentation of STEPS API, documentation relative to API_1 is still available (:doc:`API_1/API_ref`).
